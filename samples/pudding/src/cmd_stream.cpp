#include <cmd_stream.hpp>

#include <image_view.hpp>

#include "ensure.hpp"
#include "context.hpp"

using namespace pd;

render_stream::render_stream(cmd_stream& parent)
    : parent_{parent}
{
}

void cmd_stream::push_render_target(render_target const& input)
{
    rts_.emplace_back(&input);
}

cmd_stream::cmd_stream(void* cb)
    : cb_{cb}
{}

void cmd_stream::begin_render_pass()
{
    ENSURE(!rts_.empty(), "Attempting to invoke a render pass with no render targets");
    // First, determine if we've already encountered this set of render targets and retrieve it if so

    // The hash combining algorithm below is simple and order-independent.
    unique_id::hash_t rp_hash{};
    unique_id::hash_t fb_hash{};
    for (auto const* rt : rts_)
    {
        rp_hash ^= rt->uid.hash;
        fb_hash ^= rt->iv->hash();
    }

    auto rp = vk_objects.render_passes.find(rp_hash);
    if (rp == vk_objects.render_passes.end())
    {
        // Create a render pass

        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkRenderPassCreateInfo.html
        VkRenderPassCreateInfo pass_info{};
        pass_info.sType        = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
        pass_info.subpassCount = 1; // TODO

        std::vector<VkAttachmentDescription> attachments;
        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkAttachmentDescription.html
        for (auto const* rt : rts_)
        {
            VkAttachmentDescription desc{};
            desc.format = static_cast<VkFormat>(rt->format);
            desc.samples = static_cast<VkSampleCountFlagBits>(rt->sample_count_bits);

            if (rt->load)
            {
                desc.loadOp = VK_ATTACHMENT_LOAD_OP_LOAD;
                if (rt->stencil)
                {
                    desc.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_LOAD;
                }
            }
            else if (rt->clear)
            {
                desc.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
                if (rt->stencil)
                {
                    desc.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
                }
            }
            else
            {
                desc.loadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
                if (rt->stencil)
                {
                    desc.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
                }
            }

            if (rt->store)
            {
                desc.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
                if (rt->stencil)
                {
                    desc.stencilStoreOp = VK_ATTACHMENT_STORE_OP_STORE;
                }
            }
            else
            {
                desc.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
                if (rt->stencil)
                {
                    desc.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
                }
            }

            desc.initialLayout = static_cast<VkImageLayout>(rt->initial_layout);
            desc.finalLayout = static_cast<VkImageLayout>(rt->final_layout);

            attachments.emplace_back(desc);
        }

        pass_info.attachmentCount = static_cast<uint32_t>(attachments.size());
        pass_info.pAttachments = attachments.data();

        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkSubpassDescription.html
        VkSubpassDescription subpass{};
        subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
        // TODO handle input attachments
        subpass.inputAttachmentCount = 0;
        subpass.pInputAttachments    = nullptr;
        // TODO handle preserve attachments
        subpass.preserveAttachmentCount = 0;
        // TODO handle resolve attachments
        subpass.pResolveAttachments = nullptr;

        std::vector<VkAttachmentReference> color_refs;
        VkAttachmentReference ds_ref;
        uint32_t index = 0;

        bool ds_attached = false;

        // For each attachment, specify the layout they should be in *during* the render pass
        for (auto const* rt : rts_)
        {
            if (rt->depth || rt->stencil)
            {
                ENSURE(!ds_attached, "Attempted to attach more than one depth-stencil targets to a render pass");
                ds_attached = true;
                ds_ref.attachment = index++;
                ds_ref.layout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;
            }
            else
            {
                VkAttachmentReference color_ref{};
                color_ref.attachment = index++;
                color_ref.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
                color_refs.emplace_back(color_ref);
            }
        }

        subpass.pDepthStencilAttachment = ds_attached ? &ds_ref : nullptr;
        subpass.colorAttachmentCount = static_cast<uint32_t>(color_refs.size());
        subpass.pColorAttachments = color_refs.data();

        pass_info.pSubpasses = &subpass;

        // We're going to assume here that every render pass writes during the color attachment output stage and will do
        // so in multiple frames. This dependency states that before we try to write to any color attachments, wait
        // until writes prior to the start of the render pass have been flushed.
        VkSubpassDependency dep{};
        dep.srcSubpass    = VK_SUBPASS_EXTERNAL;
        dep.srcStageMask  = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dep.srcAccessMask = 0;
        dep.dstSubpass    = 0;
        dep.dstStageMask  = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        dep.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

        pass_info.dependencyCount = 1;
        pass_info.pDependencies = &dep;
        VkRenderPass render_pass;
        ENSURE(vkCreateRenderPass(vk_device, &pass_info, vk_ac, &render_pass), "Failed to create vulkan render pass");
        rp = vk_objects.render_passes.emplace(std::make_pair(rp_hash, render_pass)).first;
    }

    VkRenderPass render_pass = rp->second;

    auto fb = vk_objects.framebuffers.find(fb_hash);
    if (fb == vk_objects.framebuffers.end())
    {
        // This set of render targets has not been used before (or it was reclaimed due to misuse). Create the
        // corresponding framebuffer.

        std::vector<VkImageView> views;
        for (auto const* rt : rts_)
        {
            views.emplace_back(static_cast<VkImageView>(rt->iv->handle()));
        }

        // https://www.khronos.org/registry/vulkan/specs/1.1-extensions/man/html/VkFramebufferCreateInfo.html
        VkFramebufferCreateInfo fb_info{};
        fb_info.sType           = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
        fb_info.flags           = 0;
        fb_info.renderPass      = render_pass;
        fb_info.attachmentCount = static_cast<uint32_t>(views.size());
        fb_info.pAttachments    = views.data();
        fb_info.width           = rts_[0]->iv->width();
        fb_info.height          = rts_[0]->iv->height();
        fb_info.layers          = 1; // Rendering to texture array not supported currently

        VkFramebuffer framebuffer;
        ENSURE(vkCreateFramebuffer(vk_device, &fb_info, vk_ac, &framebuffer), "Unable to create vulkan framebuffer");
        fb = vk_objects.framebuffers.emplace(std::make_pair(fb_hash, framebuffer)).first;
    }

    VkFramebuffer framebuffer = fb->second;

    std::vector<VkClearValue> clear_values;
    for (auto* rt : rts_)
    {
        if (rt->clear)
        {
            VkClearValue cv;
            if (rt->depth || rt->stencil)
            {
                cv.depthStencil
                    = VkClearDepthStencilValue{rt->depth_stencil_clear.depth, rt->depth_stencil_clear.stencil};
            }
            else
            {
                cv.color
                    = VkClearColorValue{rt->clear_value[0], rt->clear_value[1], rt->clear_value[2], rt->clear_value[3]};
            }
            clear_values.emplace_back(cv);
        }
    }

    // We now have the necessary framebuffer and renderpass to render to. Let's start the render pass.
    VkRenderPassBeginInfo rp_info{};
    rp_info.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
    rp_info.renderArea.offset.x = 0;
    rp_info.renderArea.offset.y = 0;
    rp_info.renderArea.extent.width = rts_[0]->iv->width();
    rp_info.renderArea.extent.height = rts_[0]->iv->height();
    rp_info.clearValueCount = static_cast<uint32_t>(clear_values.size());
    rp_info.pClearValues = clear_values.data();
    rp_info.renderPass = render_pass;
    rp_info.framebuffer = framebuffer;

    vkCmdBeginRenderPass(static_cast<VkCommandBuffer>(cb_), &rp_info, VK_SUBPASS_CONTENTS_INLINE);
}

void cmd_stream::end_render_pass()
{
    vkCmdEndRenderPass(static_cast<VkCommandBuffer>(cb_));
}